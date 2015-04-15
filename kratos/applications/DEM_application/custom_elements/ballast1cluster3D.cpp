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
        
        int number_of_spheres = 34;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.0122 (in meters) was the medium diameter of the rock in GiD (Rock1_01.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.9879,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.9879;
                
        mListOfCoordinates [ 0][0] =-0.262773; mListOfCoordinates [ 0][1] = 0.019097; mListOfCoordinates [ 0][2] =-0.126412;
        mListOfCoordinates [ 1][0] = 0.036299; mListOfCoordinates [ 1][1] = 0.021229; mListOfCoordinates [ 1][2] =-0.226917;
        mListOfCoordinates [ 2][0] = 0.014782; mListOfCoordinates [ 2][1] = 0.008596; mListOfCoordinates [ 2][2] = 0.224034;
        mListOfCoordinates [ 3][0] =-0.145765; mListOfCoordinates [ 3][1] = 0.030541; mListOfCoordinates [ 3][2] =-0.287297;
        mListOfCoordinates [ 4][0] = 0.177354; mListOfCoordinates [ 4][1] = 0.006114; mListOfCoordinates [ 4][2] = 0.141173;
        mListOfCoordinates [ 5][0] =-0.238395; mListOfCoordinates [ 5][1] = 0.008703; mListOfCoordinates [ 5][2] =-0.038181;
        mListOfCoordinates [ 6][0] =-0.418637; mListOfCoordinates [ 6][1] = 0.013024; mListOfCoordinates [ 6][2] = 0.183847;
        mListOfCoordinates [ 7][0] = 0.593959; mListOfCoordinates [ 7][1] = 0.003485; mListOfCoordinates [ 7][2] = 0.014352;
        mListOfCoordinates [ 8][0] =-0.399764; mListOfCoordinates [ 8][1] = 0.022302; mListOfCoordinates [ 8][2] =-0.161183;
        mListOfCoordinates [ 9][0] =-0.349754; mListOfCoordinates [ 9][1] = 0.010286; mListOfCoordinates [ 9][2] =-0.295375;
        mListOfCoordinates [10][0] =-0.587612; mListOfCoordinates [10][1] =-0.001876; mListOfCoordinates [10][2] =-0.153889;
        mListOfCoordinates [11][0] =-0.631187; mListOfCoordinates [11][1] = 0.001949; mListOfCoordinates [11][2] = 0.064436;
        mListOfCoordinates [12][0] = 0.532218; mListOfCoordinates [12][1] = 0.011574; mListOfCoordinates [12][2] =-0.141179;
        mListOfCoordinates [13][0] =-0.056143; mListOfCoordinates [13][1] = 0.016259; mListOfCoordinates [13][2] =-0.259193;
        mListOfCoordinates [14][0] = 0.402098; mListOfCoordinates [14][1] = 0.002041; mListOfCoordinates [14][2] = 0.103020;
        mListOfCoordinates [15][0] =-0.254128; mListOfCoordinates [15][1] = 0.013865; mListOfCoordinates [15][2] = 0.248432;
        mListOfCoordinates [16][0] = 0.281280; mListOfCoordinates [16][1] = 0.031829; mListOfCoordinates [16][2] =-0.212180;
        mListOfCoordinates [17][0] =-0.393272; mListOfCoordinates [17][1] = 0.012686; mListOfCoordinates [17][2] = 0.051246;
        mListOfCoordinates [18][0] = 0.138410; mListOfCoordinates [18][1] = 0.014293; mListOfCoordinates [18][2] = 0.184005;
        mListOfCoordinates [19][0] =-0.237657; mListOfCoordinates [19][1] = 0.023515; mListOfCoordinates [19][2] =-0.249290;
        mListOfCoordinates [20][0] =-0.072404; mListOfCoordinates [20][1] =-0.001012; mListOfCoordinates [20][2] =-0.001247;
        mListOfCoordinates [21][0] = 0.175448; mListOfCoordinates [21][1] = 0.091147; mListOfCoordinates [21][2] =-0.051040;
        mListOfCoordinates [22][0] =-0.531074; mListOfCoordinates [22][1] = 0.023601; mListOfCoordinates [22][2] =-0.010766;
        mListOfCoordinates [23][0] =-0.170481; mListOfCoordinates [23][1] =-0.003579; mListOfCoordinates [23][2] = 0.068976;
        mListOfCoordinates [24][0] =-0.445309; mListOfCoordinates [24][1] = 0.021877; mListOfCoordinates [24][2] =-0.169820;
        mListOfCoordinates [25][0] =-0.150670; mListOfCoordinates [25][1] = 0.017657; mListOfCoordinates [25][2] = 0.215104;
        mListOfCoordinates [26][0] =-0.091479; mListOfCoordinates [26][1] = 0.013124; mListOfCoordinates [26][2] =-0.142943;
        mListOfCoordinates [27][0] =-0.175512; mListOfCoordinates [27][1] = 0.012355; mListOfCoordinates [27][2] = 0.217805;
        mListOfCoordinates [28][0] =-0.085496; mListOfCoordinates [28][1] =-0.002648; mListOfCoordinates [28][2] = 0.085972;
        mListOfCoordinates [29][0] = 0.229786; mListOfCoordinates [29][1] = 0.020764; mListOfCoordinates [29][2] =-0.123575;
        mListOfCoordinates [30][0] =-0.015041; mListOfCoordinates [30][1] =-0.099307; mListOfCoordinates [30][2] = 0.077238;
        mListOfCoordinates [31][0] = 0.345687; mListOfCoordinates [31][1] = 0.005533; mListOfCoordinates [31][2] = 0.033436;
        mListOfCoordinates [32][0] = 0.047961; mListOfCoordinates [32][1] = 0.022212; mListOfCoordinates [32][2] =-0.064253;
        mListOfCoordinates [33][0] = 0.351600; mListOfCoordinates [33][1] = 0.010156; mListOfCoordinates [33][2] =-0.075960;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.199183;
        mListOfRadii[ 1]= 0.218745;
        mListOfRadii[ 2]= 0.255778;
        mListOfRadii[ 3]= 0.184401;
        mListOfRadii[ 4]= 0.268336;
        mListOfRadii[ 5]= 0.209218;
        mListOfRadii[ 6]= 0.206688;
        mListOfRadii[ 7]= 0.184994;
        mListOfRadii[ 8]= 0.183692;
        mListOfRadii[ 9]= 0.194576;
        mListOfRadii[10]= 0.166181;
        mListOfRadii[11]= 0.165078;
        mListOfRadii[12]= 0.218567;
        mListOfRadii[13]= 0.181078;
        mListOfRadii[14]= 0.212356;
        mListOfRadii[15]= 0.236673;
        mListOfRadii[16]= 0.142964;
        mListOfRadii[17]= 0.192695;
        mListOfRadii[18]= 0.197730;
        mListOfRadii[19]= 0.191680;
        mListOfRadii[20]= 0.202352;
        mListOfRadii[21]= 0.100967;
        mListOfRadii[22]= 0.169001;
        mListOfRadii[23]= 0.202381;
        mListOfRadii[24]= 0.200122;
        mListOfRadii[25]= 0.200207;
        mListOfRadii[26]= 0.212194;
        mListOfRadii[27]= 0.194278;
        mListOfRadii[28]= 0.192909;
        mListOfRadii[29]= 0.169839;
        mListOfRadii[30]= 0.103225;
        mListOfRadii[31]= 0.195763;
        mListOfRadii[32]= 0.207302;
        mListOfRadii[33]= 0.191825;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.386097 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.059503047;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.181064606;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.145589529; 
  
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

