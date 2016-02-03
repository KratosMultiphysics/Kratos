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
        
        // 0.85 (in meters) was the medium diameter of the rock in GiD (Rock1_01.gid file)
        // so we should multiply every size that follow by the inverse of that number, 1.176,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 1.176;

        mListOfCoordinates[ 0][0]= 0.128142; mListOfCoordinates[ 0][1]=-0.014079; mListOfCoordinates[ 0][2]= 0.121829;
        mListOfCoordinates[ 1][0]=-0.071335; mListOfCoordinates[ 1][1]=-0.000395; mListOfCoordinates[ 1][2]= 0.266752;
        mListOfCoordinates[ 2][0]=-0.234856; mListOfCoordinates[ 2][1]= 0.012576; mListOfCoordinates[ 2][2]= 0.282531;
        mListOfCoordinates[ 3][0]=-0.441824; mListOfCoordinates[ 3][1]= 0.003334; mListOfCoordinates[ 3][2]= 0.066589;
        mListOfCoordinates[ 4][0]=-0.387133; mListOfCoordinates[ 4][1]= 0.014883; mListOfCoordinates[ 4][2]= 0.213451;
        mListOfCoordinates[ 5][0]= 0.343852; mListOfCoordinates[ 5][1]= 0.017795; mListOfCoordinates[ 5][2]=-0.219797;
        mListOfCoordinates[ 6][0]= 0.065101; mListOfCoordinates[ 6][1]=-0.004967; mListOfCoordinates[ 6][2]= 0.229315;
        mListOfCoordinates[ 7][0]=-0.231766; mListOfCoordinates[ 7][1]= 0.007212; mListOfCoordinates[ 7][2]=-0.302533;
        mListOfCoordinates[ 8][0]=-0.651388; mListOfCoordinates[ 8][1]=-0.028065; mListOfCoordinates[ 8][2]= 0.049515;
        mListOfCoordinates[ 9][0]= 0.029735; mListOfCoordinates[ 9][1]= 0.012875; mListOfCoordinates[ 9][2]=-0.250808;
        mListOfCoordinates[10][0]= 0.087384; mListOfCoordinates[10][1]= 0.000230; mListOfCoordinates[10][2]= 0.008210;
        mListOfCoordinates[11][0]= 0.550818; mListOfCoordinates[11][1]=-0.046971; mListOfCoordinates[11][2]= 0.102865;
        mListOfCoordinates[12][0]= 0.161167; mListOfCoordinates[12][1]= 0.013206; mListOfCoordinates[12][2]=-0.204964;
        mListOfCoordinates[13][0]= 0.056542; mListOfCoordinates[13][1]= 0.030440; mListOfCoordinates[13][2]= 0.312696;
        mListOfCoordinates[14][0]= 0.337199; mListOfCoordinates[14][1]=-0.018771; mListOfCoordinates[14][2]= 0.144421;
        mListOfCoordinates[15][0]=-0.093664; mListOfCoordinates[15][1]= 0.013248; mListOfCoordinates[15][2]=-0.278272;
        mListOfCoordinates[16][0]= 0.215115; mListOfCoordinates[16][1]= 0.002262; mListOfCoordinates[16][2]=-0.128458;
        mListOfCoordinates[17][0]= 0.261385; mListOfCoordinates[17][1]= 0.002196; mListOfCoordinates[17][2]=-0.042657;
        mListOfCoordinates[18][0]=-0.623520; mListOfCoordinates[18][1]=-0.027326; mListOfCoordinates[18][2]=-0.018428;
        mListOfCoordinates[19][0]=-0.601415; mListOfCoordinates[19][1]=-0.011521; mListOfCoordinates[19][2]= 0.119087;
        mListOfCoordinates[20][0]= 0.425699; mListOfCoordinates[20][1]=-0.007240; mListOfCoordinates[20][2]=-0.056806;
        mListOfCoordinates[21][0]= 0.240373; mListOfCoordinates[21][1]=-0.011582; mListOfCoordinates[21][2]= 0.164892;
        mListOfCoordinates[22][0]= 0.422456; mListOfCoordinates[22][1]=-0.015524; mListOfCoordinates[22][2]= 0.116793;
        mListOfCoordinates[23][0]= 0.184563; mListOfCoordinates[23][1]= 0.012585; mListOfCoordinates[23][2]= 0.212550;
        mListOfCoordinates[24][0]=-0.273628; mListOfCoordinates[24][1]= 0.007811; mListOfCoordinates[24][2]=-0.172399;
        mListOfCoordinates[25][0]= 0.458972; mListOfCoordinates[25][1]= 0.002935; mListOfCoordinates[25][2]=-0.128836;
        mListOfCoordinates[26][0]= 0.733911; mListOfCoordinates[26][1]=-0.001555; mListOfCoordinates[26][2]=-0.014908;
        mListOfCoordinates[27][0]=-0.523123; mListOfCoordinates[27][1]= 0.006525; mListOfCoordinates[27][2]=-0.060181;
        mListOfCoordinates[28][0]=-0.432792; mListOfCoordinates[28][1]= 0.002716; mListOfCoordinates[28][2]=-0.187953;
        mListOfCoordinates[29][0]=-0.508217; mListOfCoordinates[29][1]= 0.008538; mListOfCoordinates[29][2]= 0.184822;
        mListOfCoordinates[30][0]=-0.090562; mListOfCoordinates[30][1]=-0.018008; mListOfCoordinates[30][2]= 0.092168;
        mListOfCoordinates[31][0]=-0.332869; mListOfCoordinates[31][1]=-0.002437; mListOfCoordinates[31][2]=-0.257721;
        mListOfCoordinates[32][0]=-0.065437; mListOfCoordinates[32][1]= 0.000269; mListOfCoordinates[32][2]=-0.156860;
        mListOfCoordinates[33][0]=-0.381706; mListOfCoordinates[33][1]= 0.115076; mListOfCoordinates[33][2]= 0.322290;
        mListOfCoordinates[34][0]=-0.283287; mListOfCoordinates[34][1]=-0.010103; mListOfCoordinates[34][2]= 0.141593;
        mListOfCoordinates[35][0]= 0.572414; mListOfCoordinates[35][1]=-0.004599; mListOfCoordinates[35][2]=-0.029802;
        mListOfCoordinates[36][0]=-0.310949; mListOfCoordinates[36][1]= 0.002259; mListOfCoordinates[36][2]=-0.042018;
        mListOfCoordinates[37][0]=-0.211158; mListOfCoordinates[37][1]= 0.000533; mListOfCoordinates[37][2]=-0.080139;
        mListOfCoordinates[38][0]= 0.102414; mListOfCoordinates[38][1]= 0.004664; mListOfCoordinates[38][2]=-0.136588;
        mListOfCoordinates[39][0]= 0.634670; mListOfCoordinates[39][1]=-0.009240; mListOfCoordinates[39][2]=-0.160353;
        mListOfCoordinates[40][0]= 0.526544; mListOfCoordinates[40][1]=-0.004369; mListOfCoordinates[40][2]=-0.228498;
        mListOfCoordinates[41][0]=-0.617729; mListOfCoordinates[41][1]=-0.008737; mListOfCoordinates[41][2]= 0.223700;
        mListOfCoordinates[42][0]=-0.610101; mListOfCoordinates[42][1]=-0.009665; mListOfCoordinates[42][2]=-0.176391;
        mListOfCoordinates[43][0]= 0.706573; mListOfCoordinates[43][1]=-0.033734; mListOfCoordinates[43][2]= 0.086968;
        mListOfCoordinates[44][0]= 0.733416; mListOfCoordinates[44][1]=-0.005259; mListOfCoordinates[44][2]=-0.127669;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.214246;
        mListOfRadii[ 1]= 0.203853;
        mListOfRadii[ 2]= 0.194311;
        mListOfRadii[ 3]= 0.162074;
        mListOfRadii[ 4]= 0.157631;
        mListOfRadii[ 5]= 0.141653;
        mListOfRadii[ 6]= 0.205914;
        mListOfRadii[ 7]= 0.170504;
        mListOfRadii[ 8]= 0.102614;
        mListOfRadii[ 9]= 0.185937;
        mListOfRadii[10]= 0.210523;
        mListOfRadii[11]= 0.159393;
        mListOfRadii[12]= 0.179317;
        mListOfRadii[13]= 0.158750;
        mListOfRadii[14]= 0.212013;
        mListOfRadii[15]= 0.189565;
        mListOfRadii[16]= 0.164477;
        mListOfRadii[17]= 0.200691;
        mListOfRadii[18]= 0.088999;
        mListOfRadii[19]= 0.112286;
        mListOfRadii[20]= 0.160873;
        mListOfRadii[21]= 0.216086;
        mListOfRadii[22]= 0.208243;
        mListOfRadii[23]= 0.204754;
        mListOfRadii[24]= 0.189953;
        mListOfRadii[25]= 0.127255;
        mListOfRadii[26]= 0.081570;
        mListOfRadii[27]= 0.167172;
        mListOfRadii[28]= 0.187151;
        mListOfRadii[29]= 0.120372;
        mListOfRadii[30]= 0.245024;
        mListOfRadii[31]= 0.170602;
        mListOfRadii[32]= 0.223282;
        mListOfRadii[33]= 0.040010;
        mListOfRadii[34]= 0.186108;
        mListOfRadii[35]= 0.154576;
        mListOfRadii[36]= 0.197468;
        mListOfRadii[37]= 0.195567;
        mListOfRadii[38]= 0.183652;
        mListOfRadii[39]= 0.091947;
        mListOfRadii[40]= 0.098885;
        mListOfRadii[41]= 0.022824;
        mListOfRadii[42]= 0.100439;
        mListOfRadii[43]= 0.036420;
        mListOfRadii[44]= 0.069758;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.267371 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.048078;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.112796;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.144497;
    }
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast1Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast1Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast1Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

