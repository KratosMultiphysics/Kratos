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
#include "ballast5cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast5Cluster3D::Ballast5Cluster3D() : Cluster3D() {}
            
      
    Ballast5Cluster3D::Ballast5Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast5Cluster3D::Ballast5Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast5Cluster3D::Ballast5Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast5Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast5Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast5Cluster3D::~Ballast5Cluster3D() {}
      
    
    void Ballast5Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 45;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.1319 (in meters) was the medium diameter of the rock in GiD (grain_0451_geom.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.8835,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.8835;

        mListOfCoordinates[ 0][0]= 0.078166; mListOfCoordinates[ 0][1]= 0.097599; mListOfCoordinates[ 0][2]= 0.287078;
        mListOfCoordinates[ 1][0]=-0.037730; mListOfCoordinates[ 1][1]= 0.036351; mListOfCoordinates[ 1][2]= 0.176444;
        mListOfCoordinates[ 2][0]=-0.379995; mListOfCoordinates[ 2][1]=-0.056796; mListOfCoordinates[ 2][2]= 0.194826;
        mListOfCoordinates[ 3][0]= 0.018577; mListOfCoordinates[ 3][1]= 0.130706; mListOfCoordinates[ 3][2]= 0.463459;
        mListOfCoordinates[ 4][0]= 0.340181; mListOfCoordinates[ 4][1]=-0.088377; mListOfCoordinates[ 4][2]= 0.071023;
        mListOfCoordinates[ 5][0]=-0.093895; mListOfCoordinates[ 5][1]= 0.038560; mListOfCoordinates[ 5][2]= 0.460211;
        mListOfCoordinates[ 6][0]=-0.243792; mListOfCoordinates[ 6][1]= 0.078780; mListOfCoordinates[ 6][2]=-0.501359;
        mListOfCoordinates[ 7][0]=-0.248730; mListOfCoordinates[ 7][1]=-0.007089; mListOfCoordinates[ 7][2]= 0.482466;
        mListOfCoordinates[ 8][0]=-0.148606; mListOfCoordinates[ 8][1]= 0.088745; mListOfCoordinates[ 8][2]=-0.608336;
        mListOfCoordinates[ 9][0]=-0.026158; mListOfCoordinates[ 9][1]= 0.088897; mListOfCoordinates[ 9][2]=-0.470852;
        mListOfCoordinates[10][0]= 0.067457; mListOfCoordinates[10][1]= 0.070420; mListOfCoordinates[10][2]=-0.222745;
        mListOfCoordinates[11][0]= 0.044928; mListOfCoordinates[11][1]= 0.114358; mListOfCoordinates[11][2]= 0.132150;
        mListOfCoordinates[12][0]=-0.077276; mListOfCoordinates[12][1]= 0.098139; mListOfCoordinates[12][2]=-0.281283;
        mListOfCoordinates[13][0]=-0.297908; mListOfCoordinates[13][1]= 0.030517; mListOfCoordinates[13][2]=-0.069897;
        mListOfCoordinates[14][0]= 0.349751; mListOfCoordinates[14][1]=-0.097380; mListOfCoordinates[14][2]=-0.063883;
        mListOfCoordinates[15][0]=-0.203307; mListOfCoordinates[15][1]=-0.003661; mListOfCoordinates[15][2]= 0.312526;
        mListOfCoordinates[16][0]=-0.067138; mListOfCoordinates[16][1]=-0.038810; mListOfCoordinates[16][2]= 0.290697;
        mListOfCoordinates[17][0]=-0.357473; mListOfCoordinates[17][1]=-0.056505; mListOfCoordinates[17][2]= 0.421100;
        mListOfCoordinates[18][0]=-0.199730; mListOfCoordinates[18][1]= 0.030582; mListOfCoordinates[18][2]= 0.048574;
        mListOfCoordinates[19][0]=-0.214674; mListOfCoordinates[19][1]=-0.004388; mListOfCoordinates[19][2]= 0.244994;
        mListOfCoordinates[20][0]=-0.402518; mListOfCoordinates[20][1]=-0.097300; mListOfCoordinates[20][2]= 0.286447;
        mListOfCoordinates[21][0]=-0.039011; mListOfCoordinates[21][1]= 0.010501; mListOfCoordinates[21][2]=-0.611813;
        mListOfCoordinates[22][0]= 0.255684; mListOfCoordinates[22][1]= 0.077095; mListOfCoordinates[22][2]= 0.386282;
        mListOfCoordinates[23][0]= 0.302970; mListOfCoordinates[23][1]=-0.116853; mListOfCoordinates[23][2]= 0.179669;
        mListOfCoordinates[24][0]=-0.194695; mListOfCoordinates[24][1]= 0.085765; mListOfCoordinates[24][2]=-0.408285;
        mListOfCoordinates[25][0]= 0.225240; mListOfCoordinates[25][1]=-0.030871; mListOfCoordinates[25][2]= 0.283052;
        mListOfCoordinates[26][0]=-0.021730; mListOfCoordinates[26][1]= 0.079703; mListOfCoordinates[26][2]=-0.076048;
        mListOfCoordinates[27][0]=-0.338819; mListOfCoordinates[27][1]= 0.003853; mListOfCoordinates[27][2]= 0.069651;
        mListOfCoordinates[28][0]= 0.222253; mListOfCoordinates[28][1]=-0.154448; mListOfCoordinates[28][2]= 0.246453;
        mListOfCoordinates[29][0]=-0.196151; mListOfCoordinates[29][1]= 0.094630; mListOfCoordinates[29][2]=-0.295711;
        mListOfCoordinates[30][0]= 0.182386; mListOfCoordinates[30][1]=-0.031656; mListOfCoordinates[30][2]= 0.036677;
        mListOfCoordinates[31][0]= 0.370456; mListOfCoordinates[31][1]=-0.122675; mListOfCoordinates[31][2]=-0.171530;
        mListOfCoordinates[32][0]= 0.239665; mListOfCoordinates[32][1]=-0.053048; mListOfCoordinates[32][2]=-0.124781;
        mListOfCoordinates[33][0]= 0.077247; mListOfCoordinates[33][1]= 0.083528; mListOfCoordinates[33][2]=-0.374975;
        mListOfCoordinates[34][0]= 0.005977; mListOfCoordinates[34][1]=-0.159240; mListOfCoordinates[34][2]= 0.305982;
        mListOfCoordinates[35][0]= 0.082314; mListOfCoordinates[35][1]= 0.061031; mListOfCoordinates[35][2]=-0.081658;
        mListOfCoordinates[36][0]= 0.175360; mListOfCoordinates[36][1]=-0.033514; mListOfCoordinates[36][2]= 0.172601;
        mListOfCoordinates[37][0]= 0.186398; mListOfCoordinates[37][1]=-0.020377; mListOfCoordinates[37][2]=-0.412900;
        mListOfCoordinates[38][0]= 0.307012; mListOfCoordinates[38][1]=-0.059483; mListOfCoordinates[38][2]=-0.350576;
        mListOfCoordinates[39][0]= 0.045948; mListOfCoordinates[39][1]=-0.008517; mListOfCoordinates[39][2]=-0.522712;
        mListOfCoordinates[40][0]= 0.222741; mListOfCoordinates[40][1]=-0.025492; mListOfCoordinates[40][2]=-0.301808;
        mListOfCoordinates[41][0]=-0.203048; mListOfCoordinates[41][1]= 0.066359; mListOfCoordinates[41][2]=-0.164350;
        mListOfCoordinates[42][0]= 0.064476; mListOfCoordinates[42][1]=-0.126083; mListOfCoordinates[42][2]= 0.253539;
        mListOfCoordinates[43][0]=-0.216912; mListOfCoordinates[43][1]= 0.016122; mListOfCoordinates[43][2]= 0.570750;
        mListOfCoordinates[44][0]= 0.344131; mListOfCoordinates[44][1]=-0.089691; mListOfCoordinates[44][2]=-0.261133;

        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.204898;
        mListOfRadii[ 1]= 0.248121;
        mListOfRadii[ 2]= 0.129997;
        mListOfRadii[ 3]= 0.110971;
        mListOfRadii[ 4]= 0.101141;
        mListOfRadii[ 5]= 0.18327;
        mListOfRadii[ 6]= 0.076334;
        mListOfRadii[ 7]= 0.126306;
        mListOfRadii[ 8]= 0.158541;
        mListOfRadii[ 9]= 0.203504;
        mListOfRadii[10]= 0.211509;
        mListOfRadii[11]= 0.207299;
        mListOfRadii[12]= 0.154789;
        mListOfRadii[13]= 0.084444;
        mListOfRadii[14]= 0.090622;
        mListOfRadii[15]= 0.205679;
        mListOfRadii[16]= 0.219614;
        mListOfRadii[17]= 0.096091;
        mListOfRadii[18]= 0.156978;
        mListOfRadii[19]= 0.233531;
        mListOfRadii[20]= 0.094059;
        mListOfRadii[21]= 0.101306;
        mListOfRadii[22]= 0.049709;
        mListOfRadii[23]= 0.116412;
        mListOfRadii[24]= 0.113252;
        mListOfRadii[25]= 0.162797;
        mListOfRadii[26]= 0.213357;
        mListOfRadii[27]= 0.124343;
        mListOfRadii[28]= 0.139851;
        mListOfRadii[29]= 0.097662;
        mListOfRadii[30]= 0.168403;
        mListOfRadii[31]= 0.077644;
        mListOfRadii[32]= 0.150284;
        mListOfRadii[33]= 0.199167;
        mListOfRadii[34]= 0.124836;
        mListOfRadii[35]= 0.209622;
        mListOfRadii[36]= 0.186548;
        mListOfRadii[37]= 0.138241;
        mListOfRadii[38]= 0.085285;
        mListOfRadii[39]= 0.13044;
        mListOfRadii[40]= 0.128303;
        mListOfRadii[41]= 0.112418;
        mListOfRadii[42]= 0.15351;
        mListOfRadii[43]= 0.114211;
        mListOfRadii[44]= 0.073588;
        
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.251605 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.045838;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.113205;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.135844;

    }
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast5Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast5Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast5Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

