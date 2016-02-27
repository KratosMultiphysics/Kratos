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
#include "ballast1cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast1Cluster3D::Ballast1Cluster3D() : Cluster3D() {}
            
      
    Ballast1Cluster3D::Ballast1Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast1Cluster3D::Ballast1Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast1Cluster3D::Ballast1Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast1Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast1Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast1Cluster3D::~Ballast1Cluster3D() {}
      
    
    void Ballast1Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 45;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 0.856789 (in meters) was the medium diameter of the rock in GiD (Rock1_01.gid file)
        // so we should multiply every size that follow by the inverse of that number, 1.167149,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 1.167149;

        mListOfCoordinates[ 0][0]=-0.409733; mListOfCoordinates[ 0][1]= 0.123811; mListOfCoordinates[ 0][2]= 0.297569;
        mListOfCoordinates[ 1][0]= 0.330120; mListOfCoordinates[ 1][1]= 0.018637; mListOfCoordinates[ 1][2]=-0.201518;
        mListOfCoordinates[ 2][0]=-0.604839; mListOfCoordinates[ 2][1]=-0.014783; mListOfCoordinates[ 2][2]= 0.044262;
        mListOfCoordinates[ 3][0]=-0.399951; mListOfCoordinates[ 3][1]= 0.013809; mListOfCoordinates[ 3][2]= 0.189746;
        mListOfCoordinates[ 4][0]= 0.062246; mListOfCoordinates[ 4][1]=-0.004036; mListOfCoordinates[ 4][2]= 0.209702;
        mListOfCoordinates[ 5][0]= 0.032805; mListOfCoordinates[ 5][1]=-0.014702; mListOfCoordinates[ 5][2]= 0.066936;
        mListOfCoordinates[ 6][0]=-0.062094; mListOfCoordinates[ 6][1]= 0.012667; mListOfCoordinates[ 6][2]=-0.284208;
        mListOfCoordinates[ 7][0]=-0.167470; mListOfCoordinates[ 7][1]= 0.017194; mListOfCoordinates[ 7][2]=-0.296667;
        mListOfCoordinates[ 8][0]= 0.046779; mListOfCoordinates[ 8][1]= 0.008670; mListOfCoordinates[ 8][2]=-0.250713;
        mListOfCoordinates[ 9][0]=-0.241507; mListOfCoordinates[ 9][1]= 0.005545; mListOfCoordinates[ 9][2]=-0.325262;
        mListOfCoordinates[10][0]= 0.574812; mListOfCoordinates[10][1]=-0.020080; mListOfCoordinates[10][2]= 0.061735;
        mListOfCoordinates[11][0]= 0.676110; mListOfCoordinates[11][1]= 0.002623; mListOfCoordinates[11][2]=-0.031214;
        mListOfCoordinates[12][0]= 0.461822; mListOfCoordinates[12][1]= 0.002886; mListOfCoordinates[12][2]=-0.125331;
        mListOfCoordinates[13][0]= 0.118720; mListOfCoordinates[13][1]= 0.010866; mListOfCoordinates[13][2]=-0.087631;
        mListOfCoordinates[14][0]= 0.294502; mListOfCoordinates[14][1]= 0.002858; mListOfCoordinates[14][2]=-0.084549;
        mListOfCoordinates[15][0]=-0.402572; mListOfCoordinates[15][1]=-0.001231; mListOfCoordinates[15][2]=-0.221673;
        mListOfCoordinates[16][0]=-0.249441; mListOfCoordinates[16][1]= 0.026585; mListOfCoordinates[16][2]= 0.294170;
        mListOfCoordinates[17][0]=-0.376636; mListOfCoordinates[17][1]=-0.002807; mListOfCoordinates[17][2]= 0.123651;
        mListOfCoordinates[18][0]= 0.446462; mListOfCoordinates[18][1]=-0.009578; mListOfCoordinates[18][2]=-0.009702;
        mListOfCoordinates[19][0]= 0.197851; mListOfCoordinates[19][1]=-0.001783; mListOfCoordinates[19][2]= 0.173208;
        mListOfCoordinates[20][0]=-0.549442; mListOfCoordinates[20][1]= 0.002232; mListOfCoordinates[20][2]= 0.177873;
        mListOfCoordinates[21][0]=-0.567578; mListOfCoordinates[21][1]=-0.011804; mListOfCoordinates[21][2]=-0.170920;
        mListOfCoordinates[22][0]=-0.471918; mListOfCoordinates[22][1]= 0.011566; mListOfCoordinates[22][2]=-0.110496;
        mListOfCoordinates[23][0]=-0.181913; mListOfCoordinates[23][1]=-0.007914; mListOfCoordinates[23][2]= 0.000652;
        mListOfCoordinates[24][0]=-0.257889; mListOfCoordinates[24][1]= 0.010313; mListOfCoordinates[24][2]=-0.197178;
        mListOfCoordinates[25][0]=-0.222339; mListOfCoordinates[25][1]=-0.009515; mListOfCoordinates[25][2]= 0.169483;
        mListOfCoordinates[26][0]=-0.057411; mListOfCoordinates[26][1]= 0.004237; mListOfCoordinates[26][2]= 0.252195;
        mListOfCoordinates[27][0]=-0.413748; mListOfCoordinates[27][1]= 0.006967; mListOfCoordinates[27][2]= 0.026494;
        mListOfCoordinates[28][0]= 0.343975; mListOfCoordinates[28][1]=-0.010430; mListOfCoordinates[28][2]= 0.129125;
        mListOfCoordinates[29][0]= 0.200368; mListOfCoordinates[29][1]=-0.003267; mListOfCoordinates[29][2]= 0.051256;
        mListOfCoordinates[30][0]= 0.749982; mListOfCoordinates[30][1]=-0.005989; mListOfCoordinates[30][2]=-0.156702;
        mListOfCoordinates[31][0]=-0.069275; mListOfCoordinates[31][1]=-0.011946; mListOfCoordinates[31][2]= 0.112387;
        mListOfCoordinates[32][0]= 0.172462; mListOfCoordinates[32][1]= 0.016683; mListOfCoordinates[32][2]=-0.217875;
        mListOfCoordinates[33][0]=-0.318167; mListOfCoordinates[33][1]= 0.006446; mListOfCoordinates[33][2]=-0.077170;
        mListOfCoordinates[34][0]=-0.343545; mListOfCoordinates[34][1]= 0.001005; mListOfCoordinates[34][2]=-0.247658;
        mListOfCoordinates[35][0]= 0.084700; mListOfCoordinates[35][1]= 0.036717; mListOfCoordinates[35][2]= 0.298939;
        mListOfCoordinates[36][0]= 0.464462; mListOfCoordinates[36][1]= 0.008158; mListOfCoordinates[36][2]=-0.264878;
        mListOfCoordinates[37][0]= 0.453043; mListOfCoordinates[37][1]=-0.022546; mListOfCoordinates[37][2]= 0.109390;
        mListOfCoordinates[38][0]=-0.083184; mListOfCoordinates[38][1]= 0.000899; mListOfCoordinates[38][2]=-0.123891;
        mListOfCoordinates[39][0]= 0.685940; mListOfCoordinates[39][1]= 0.004548; mListOfCoordinates[39][2]=-0.070893;
        mListOfCoordinates[40][0]= 0.617578; mListOfCoordinates[40][1]= 0.002957; mListOfCoordinates[40][2]=-0.102599;
        mListOfCoordinates[41][0]= 0.754596; mListOfCoordinates[41][1]=-0.009442; mListOfCoordinates[41][2]= 0.011713;
        mListOfCoordinates[42][0]= 0.697203; mListOfCoordinates[42][1]=-0.010210; mListOfCoordinates[42][2]=-0.199289;
        mListOfCoordinates[43][0]= 0.771453; mListOfCoordinates[43][1]=-0.000072; mListOfCoordinates[43][2]=-0.054340;
        mListOfCoordinates[44][0]= 0.571263; mListOfCoordinates[44][1]=-0.008049; mListOfCoordinates[44][2]=-0.224916;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.020672;
        mListOfRadii[ 1]= 0.155162;
        mListOfRadii[ 2]= 0.141522;
        mListOfRadii[ 3]= 0.155502;
        mListOfRadii[ 4]= 0.204192;
        mListOfRadii[ 5]= 0.223678;
        mListOfRadii[ 6]= 0.192530;
        mListOfRadii[ 7]= 0.172583;
        mListOfRadii[ 8]= 0.189502;
        mListOfRadii[ 9]= 0.159412;
        mListOfRadii[10]= 0.178333;
        mListOfRadii[11]= 0.076550;
        mListOfRadii[12]= 0.147031;
        mListOfRadii[13]= 0.206146;
        mListOfRadii[14]= 0.190871;
        mListOfRadii[15]= 0.171714;
        mListOfRadii[16]= 0.181904;
        mListOfRadii[17]= 0.146873;
        mListOfRadii[18]= 0.161730;
        mListOfRadii[19]= 0.231182;
        mListOfRadii[20]= 0.119456;
        mListOfRadii[21]= 0.137786;
        mListOfRadii[22]= 0.185702;
        mListOfRadii[23]= 0.220965;
        mListOfRadii[24]= 0.200797;
        mListOfRadii[25]= 0.206410;
        mListOfRadii[26]= 0.219153;
        mListOfRadii[27]= 0.175220;
        mListOfRadii[28]= 0.220325;
        mListOfRadii[29]= 0.198387;
        mListOfRadii[30]= 0.050396;
        mListOfRadii[31]= 0.196537;
        mListOfRadii[32]= 0.179728;
        mListOfRadii[33]= 0.178625;
        mListOfRadii[34]= 0.174406;
        mListOfRadii[35]= 0.157729;
        mListOfRadii[36]= 0.068093;
        mListOfRadii[37]= 0.198215;
        mListOfRadii[38]= 0.222955;
        mListOfRadii[39]= 0.073831;
        mListOfRadii[40]= 0.115712;
        mListOfRadii[41]= 0.039945;
        mListOfRadii[42]= 0.031858;
        mListOfRadii[43]= 0.065943;
        mListOfRadii[44]= 0.103893;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.329321 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.080199;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.143777;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.167324;
    }
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast1Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Ballast1Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Ballast1Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

